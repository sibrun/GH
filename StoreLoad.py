import os
import pickle


class NotBuiltError(RuntimeError):
    pass


class RefError(RuntimeError):
    pass


class FileNotExistingError(RuntimeError):
    pass


def generate_path(path):
    directory = os.path.dirname(path)
    if not os.path.exists(directory):
        os.makedirs(directory)


def store_string_list(L, path):
    generate_path(path)
    with open(path,'w') as f:
        for x in L:
            f.write(x + '\n')


def load_string_list(path):
    if not os.path.exists(path):
        raise FileNotExistingError("Cannot load from %s: The file does not exist" % str(path))
    with open(path, 'r') as f:
        return f.read().splitlines()


def load_line(path):
    if not os.path.exists(path):
        raise FileNotExistingError("Cannot load from %s: The file does not exist" % str(path))
    with open(path, 'r') as f:
        return f.readline()


def store_line(S, path):
    generate_path(path)
    with open(path,'w') as f:
        f.write(S + '\n')


def store_list_of_header_lists(LHL, path):
    generate_path(path)
    with open(path, 'w') as f:
        for HL in LHL:
            (H,L) = HL
            if H is not None:
                f.write(H + '\n')
            for x in L:
                f.write(x + '\n')


def pickle_store(Ob, path):
    generate_path(path)
    with open(path,'wb') as f:
        pickle.dump(Ob, f)


def pickle_load(path):
    if not os.path.exists(path):
        raise FileNotExistingError("Cannot load from %s: The file does not exist" % str(path))
    with open(path, 'rb') as f:
        return pickle.load(f)