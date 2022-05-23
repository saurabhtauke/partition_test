import matplotlib.pyplot as plt
import numpy as np
import bisect
import time
from numba import cuda, float32
import math
import os
from scipy.io import savemat

cuda.select_device(0)
print(f'Currently using the {cuda.get_current_device()}')

def load_csv_file(target_filename):
    data = np.loadtxt(target_filename, delimiter=',', skiprows=2, usecols=[3, 8])#[:157000,:]
    print(f'Loaded csv file with shape {data.shape}')
    dtimes = data[:, 1].astype(np.uint16)
    timetags = data[:, 0].astype(np.uint64)
    return dtimes, timetags

def slow_cuda_for_one_bigdelay(target_filename = 'data/PMA_PEG_A5001_x27.80µm_y31.00µm_z40..csv',
    window=2000000,
    bigdelay = 3,
    pair_matrix_size=False,
    TPB=8,
    do_plot=True):

    t0 = time.time()
    dtimes, timetags = load_csv_file(target_filename)
    mintimetag = timetags[0]
    maxtimetag = timetags[-1]
    print(f'max timetag is {np.max(timetags)}')
    print(f'This means number of windows = {(maxtimetag - mintimetag) // window}')

    if not pair_matrix_size:
        pair_matrix_size = np.max(dtimes) + 1
        print(f'Pair matrix auto size: {pair_matrix_size}')

    # split the data into batches.
    # Drop the last batch
    # TODO: This may be sped up by assuming the linear approximation between the index and time and using np.reshape instead
    batch_start_timetags = np.arange(mintimetag, maxtimetag, window)
    nbatches = batch_start_timetags.shape[0] - 1
    batch_start_indices = [bisect.bisect_left(timetags, batch_start_timetag) for batch_start_timetag in batch_start_timetags]

    # sort the dtimes within each batch -- as a prayer to runtime CPU optimization
    for batch_id in range(nbatches):
        dtimes[batch_start_indices[batch_id]: batch_start_indices[batch_id + 1]] = \
            np.sort(dtimes[batch_start_indices[batch_id]: batch_start_indices[batch_id + 1]])

    # bigdelay_in_windows = int(round(big_delay_overall / window))

    pair_matrix = np.zeros(shape=(pair_matrix_size, pair_matrix_size), dtype=np.int16)
    batch_start_indices = np.array(batch_start_indices)
    threadsperblock = (TPB, TPB)
    blockspergrid_x = math.ceil(pair_matrix.shape[0] / threadsperblock[0])
    blockspergrid_y = math.ceil(pair_matrix.shape[1] / threadsperblock[1])
    blockspergrid = (blockspergrid_x, blockspergrid_y)

    t0 = time.time()
    pair_matrix_d = cuda.to_device(pair_matrix)
    dtimes_d = cuda.to_device(dtimes)
    batch_start_indices_d = cuda.to_device(batch_start_indices)

    ndelays = 1
    print(bigdelay)
    slow_cuda_thread_for_pairmatrix[blockspergrid, threadsperblock](pair_matrix_d, dtimes_d, batch_start_indices_d,
                                                                    bigdelay)
    # np.save(f'output/{bigdelay:08d}.npy', pair_matrix_d.copy_to_host())
    cuda.synchronize()
    pair_matrix = pair_matrix_d.copy_to_host()
    print(f'Cuda computing time: {(time.time() - t0):.2f} seconds')
    # plt.imshow()
    cuda.close()

    if do_plot:
        plt.imshow(pair_matrix)
        plt.colorbar()
        plt.show()

    return pair_matrix

@cuda.jit('void(uint16[:, :], uint16[:], uint32[:], uint32)')
def slow_cuda_thread_for_pairmatrix(pair_matrix, dtimes, batch_start_indices, bigdelay):
    nbatches = batch_start_indices.shape[0] - 1
    x, y = cuda.grid(2)
    if x < pair_matrix.shape[0] and y < pair_matrix.shape[1]:
        pair_matrix[x, y] = 0
        # an_array[x, y] = nbatches
        for first_batch_id in range(nbatches - bigdelay):
            first_batch_dtimes = dtimes[batch_start_indices[first_batch_id]: batch_start_indices[first_batch_id+1]]
            first_batch_count = 0
            for i in first_batch_dtimes:
                if i == x:
                    first_batch_count += 1
                elif i > x:
                    break
            if first_batch_count > 0:
                second_batch_dtimes = dtimes[batch_start_indices[first_batch_id + bigdelay]:
                                             batch_start_indices[first_batch_id + bigdelay + 1]]
                overall_count = 0
                for j in second_batch_dtimes:
                    if j == y:
                        overall_count += first_batch_count
                    elif j > y:
                        break
                if overall_count > 0:
                    pair_matrix[x, y] += overall_count

def faster_cuda_with_precounting(target_filename = 'data/PMA_PEG_A5001_x27.80µm_y31.00µm_z40..csv',
    window=2000000,
    bigdelays = 'all',
    pair_matrix_size=False,
    TPB=8,
    do_plot=False):

    t0 = time.time()
    dtimes, timetags = load_csv_file(target_filename)
    mintimetag = timetags[0]
    maxtimetag = timetags[-1]
    print(f'Min timetag is {mintimetag}')
    print(f'Max timetag is {maxtimetag}')
    print(f'This means number of windows is {(maxtimetag - mintimetag) // window}')

    if not pair_matrix_size:
        pair_matrix_size = np.max(dtimes) + 1
        print(f'Pair matrix auto size: {pair_matrix_size}')

    if bigdelays == 'all':
        bigdelays = np.arange(int(round((maxtimetag - mintimetag) // window)))

    # split the data into batches.
    # Drop the last batch
    # TODO: This may be sped up by assuming the linear approximation between the index and time and using np.reshape instead
    batch_start_timetags = np.arange(mintimetag, maxtimetag, window)
    nbatches = batch_start_timetags.shape[0] - 1
    batch_start_indices = [bisect.bisect_left(timetags, batch_start_timetag) for batch_start_timetag in batch_start_timetags]

    # replace the list of dtimes with a matrix of shape = (pair_matrix_size, nbatches) filles with counts
    # of photons having a given dtime within a given batch.
    t0batchmatrix = time.time()
    batchmatrix = np.zeros(shape=(pair_matrix_size, nbatches), dtype=np.uint16)
    for batch_id in range(nbatches):
        np.add.at(batchmatrix[:, batch_id],
                  dtimes[batch_start_indices[batch_id]: batch_start_indices[batch_id + 1]],
                  1)
    print(f'Filled batchmatrix in {(time.time() - t0batchmatrix):.2f} seconds')

    pair_matrix = np.zeros(shape=(pair_matrix_size, pair_matrix_size), dtype=np.int16)
    batch_start_indices = np.array(batch_start_indices)
    threadsperblock = (TPB, TPB)
    blockspergrid_x = math.ceil(pair_matrix.shape[0] / threadsperblock[0])
    blockspergrid_y = math.ceil(pair_matrix.shape[1] / threadsperblock[1])
    blockspergrid = (blockspergrid_x, blockspergrid_y)

    t0 = time.time()
    pair_matrix_d = cuda.to_device(pair_matrix)
    batchmatrix_d = cuda.to_device(batchmatrix)

    pair_matrices = []
    for bigdelay in bigdelays:
        if bigdelay % 100 == 0:
            print(f'Computing bigdelay {bigdelay}')
        slow_cuda_thread_for_pairmatrix_with_precounting[blockspergrid, threadsperblock](pair_matrix_d, batchmatrix_d,
                                                                        bigdelay)
        cuda.synchronize()
        pair_matrices.append(pair_matrix_d.copy_to_host())
    pair_matrices = np.array(pair_matrices)
    print(f'Cuda computing time: {(time.time() - t0):.2f} seconds for {bigdelays.shape[0]} bigdelays')
    print(f'Pair matrices data cube shape is {pair_matrices.shape}')
    cuda.close()

    if do_plot:
        plt.imshow(pair_matrices[0])
        plt.colorbar()
        plt.show()

    return pair_matrices

@cuda.jit('void(uint16[:, :], uint16[:, :], uint32)')
def slow_cuda_thread_for_pairmatrix_with_precounting(pair_matrix, batchmatrix, bigdelay):
    nbatches = batchmatrix.shape[1]
    # pair_matrix_size = pair_matrix.shape[0]
    x, y = cuda.grid(2)
    if x < pair_matrix.shape[0] and y < pair_matrix.shape[1]:
        tmp = 0
        for first_batch_id in range(nbatches - bigdelay):
            tmp += batchmatrix[x, first_batch_id] * batchmatrix[y, first_batch_id + bigdelay]
        pair_matrix[x, y] = tmp

if __name__ == '__main__':
    # slow_for_one_bigdelay()

    # target_filename = 'data/2opercent_glyPBS_2nd_2_1_1_1.csv'
    target_filename = 'data/FCS_10percent_1_1_1_1.ptu.csv'
    bigdelay = 3

    # # TESTING SLOW GPU VERSION
    # pair_matrix = slow_cuda_for_one_bigdelay(target_filename=target_filename,
    #                            window=2000000,
    #                            bigdelay=bigdelay,
    #                            pair_matrix_size=False,
    #                            do_plot=True)
    #
    # # compare correctness of GPU and CPU versions
    # pair_matrix_cpu = np.load(f'output/{bigdelay:08d}.npy')
    # print(f'Assertion that gpu is accurate as cpu: {np.isclose(pair_matrix, pair_matrix_cpu).all()}')

    #TESTING FASTER GPU VERSION
    # bigdelays = np.arange(400)
    bigdelays = 'all'
    pair_matrices = faster_cuda_with_precounting(target_filename=target_filename,
                               window=2000000,
                               bigdelays=bigdelays,
                               pair_matrix_size=False,
                               do_plot=False)

    # Save to matlab files in chunks of chunk_size, where chunk_size is number of bigdelays
    # TODO: Do parallelized version of this to speed it up, one chunk per thread
    chunk_size = 200
    do_compression = True # makes files 3 times smaller but it takes 5 times more time
    matlab_filename_prefix = os.path.splitext(target_filename)[0]
    t0 = time.time()
    chunk_indices = np.arange(0, pair_matrices.shape[0], step=chunk_size)
    for chunk_id, chunk_start_index in enumerate(chunk_indices):
        output_filename = f'{matlab_filename_prefix}_pair_matrices_chunk_{chunk_id:05d}.mat'
        print(f'Saving chunk {chunk_id} to matlab file {output_filename}')
        savemat(output_filename,
                {"number_of_bigdelays": pair_matrices.shape[0],
                 "number_of_chunks": chunk_indices.shape[0],
                 "this_chunk_id":chunk_id,
                 "pair_matrices": pair_matrices[chunk_start_index:chunk_start_index + chunk_size]},
                do_compression=do_compression)
    print(f'Saved to matlab file in {(time.time() - t0):.2f} seconds.')

    # # # compare correctness of GPU and CPU versions
    # pair_matrices = faster_cuda_with_precounting(target_filename=target_filename,
    #                            window=2000000,
    #                            bigdelays=np.array([3]),
    #                            pair_matrix_size=False,
    #                            do_plot=True)
    # pair_matrix_cpu = np.load(f'output/{3:08d}.npy')
    # print(f'Assertion that gpu is accurate as cpu: {np.isclose(pair_matrices[0], pair_matrix_cpu).all()}')

