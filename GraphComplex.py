from abc import ABCMeta, abstractmethod
import Display


class GraphComplex(object):
    __metaclass__ = ABCMeta

    def __init__(self, vector_space, operator_list):
        self.vector_space = vector_space
        self.operator_list = operator_list

    def get_vector_space(self):
        return self.vector_space

    def get_differential_list(self):
        return self.operator_list

    def build_basis(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        self.vector_space.build_basis(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs,
                                      progress_bar=progress_bar)

    def build_matrix(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        for dif in self.operator_list:
            dif.build_matrix(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs,
                                       progress_bar=progress_bar)

    def square_zero_test(self):
        for dif in self.operator_list:
            dif.square_zero_test()

    def compute_rank(self, exact=False, n_primes=1, estimate=True, ignore_existing_files=True, n_jobs=1):
        for dif in self.operator_list:
            dif.compute_rank(exact=exact, n_primes=n_primes, estimate=estimate,
                                       ignore_existing_files=ignore_existing_files, n_jobs=n_jobs)

    def plot_cohomology_dim(self, dif_idx):
        dif = self.operator_list[dif_idx]
        ordered_param_range_dict = self.vector_space.get_ordered_param_range_dict()
        dif.plot_cohomology_dim(ordered_param_range_dict)

    def commute(self, op_list1, op_list2, anti_commute=False):
        succ = []  # holds pairs for which test was successful
        fail = []  # failed pairs
        triv = []  # pairs for which test trivially succeeded because at least one operator is the empty matrix
        inc = []  # pairs for which operator matrices are missing

        for op1a in op_list1:
            for op2a in op_list2:
                if op1a.get_domain() == op2a.get_domain():
                    for op1b in op_list1:
                        if op1b.get_domain() == op2a.get_target():
                            for op2b in op_list2:
                                if op2b.get_domain() == op1a.get_target() and op1b.get_target() == op2b.get_target():
                                    p = (op1a, op1b, op2a, op2b)

                                    if not((op1a.is_valid() and op2b.is_valid()) or (op2a.is_valid() and op2b.is_valid())):
                                        triv.append(p)
                                        continue


                          '''  elseif(is_valid_op(op1a) & & is_valid_op(op2b)) & & !(
                                        is_valid_op(op1b) & & is_valid_op(op2a))
                            D = []
                            DD = []
                            try
                                D = load_matrix(op1a)
                                DD = load_matrix(op2b)
                            catch
                            # println("cannot load")
                            push!(inc, p)
                            continue
                        end
                        if D == [] | | DD == []
                            push!(triv, p)
                        else
                            if mynorm(D * DD) < 1e-10
                                push!(succ, p)
                            else
                                push!(fail, p)
                            end
                        end
                    elseif !(is_valid_op(op1a) & & is_valid_op(op2b)) & & (is_valid_op(op1b) & & is_valid_op(op2a))
                    D = []
                    DD = []
                    try
                        D = load_matrix(op2a)
                        DD = load_matrix(op1b)
                    catch
                    # println("cannot load")
                    push!(inc, p)
                    continue
                end
                if D == [] | | DD == []
                    push!(triv, p)
                else
                    if mynorm(D * DD) < 1e-10
                        push!(succ, p)
                    else
                        push!(fail, p)
                    end
                end
            else
                D1a = []
                D1b = []
                D2a = []
                D2b = []
                try
                    D1a = load_matrix(op1a)
                    D1b = load_matrix(op1b)
                    D2a = load_matrix(op2a)
                    D2b = load_matrix(op2b)
                catch
                # println("cannot load")
                push!(inc, p)
                continue
            end
            if (D1a == [] | | D2b == []) & & (D2a == [] | | D1b == [])
                push!(triv, p)
            elseif !(D1a == [] | | D2b == []) & & (D2a == [] | | D1b == [])
            if mynorm(D1a * D2b) < 1e-10
                push!(succ, p)
            else
                push!(fail, p)
            end
        elseif(D1a == [] | | D2b == []) & & !(D2a == [] | | D1b == [])
        if mynorm(D2a * D1b) < 1e-10
            push!(succ, p)
        else
            push!(fail, p)
        end
    else
        if mynorm(D1a * D2b + (antiCommute?1:-1) * D2a * D1b) < 1e-10
            push!(succ, p)
        else
            push!(fail, p)'''