import HairyGraphBiComplex as HGBC
import multiprocessing as mp

if __name__ == "__main__":

    gc = HGBC.HairyBiGC(range(6, 12), 0, True, False)

    gc.build_basis(n_jobs=1, info_tracker=True)


    '''a = [[1,2,3,4,5],[2,3,4,5]]

    def f(x):
        return x * x

    def g(aa):
        for x in aa:
            p = mp.Process(target=f,args=(x,))
        #result=[pool.apply_async(f, args=(x,)) for x in aa]
        #pool.close()
        #pool.join()
        p.start()
        p.join()
        print(aa)


    pool = mp.Pool(2)
    [pool.apply_async(g, args=(x,)) for x in a]
    pool.close()
    pool.join()'''