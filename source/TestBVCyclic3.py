import BVCyclic
import OrdinaryGraphComplex
import GraphVectorSpace
import GraphOperator

nr_jobs = 10
print(f"Building GOneVS3 bases using {nr_jobs} jobs ...")
sumvs = GraphVectorSpace.SumVectorSpace( [ BVCyclic.GOneVS3V(v, l) for v in range(0,20) for l in range(0,9) ] )

sumvs.build_basis(n_jobs=nr_jobs)

print("Finished computing bases.")

op_list = [ BVCyclic.ReconnectEdgesGO3V.generate_operator(v, l) for v in range(0,20) for l in range(0,9) ]

allop = GraphOperator.OperatorMatrixCollection(sumvs, op_list)

allop.build_matrix(n_jobs=nr_jobs)

print("Finished computing matrices.")

print("computing ranks")
allop.compute_rank(linbox="rational", n_jobs=nr_jobs)
print("Finished")