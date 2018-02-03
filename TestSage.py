from abc import ABCMeta, abstractmethod

class Fahrzeug:
    def __init__(self, w, number):
        self.wheels = w
        self.number = number


class Auto(Fahrzeug):
    def __init__(self, n):
        super(Auto, self).__init__(4, n)


class Fuhrpark:
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_fz(self, n):
        pass

    def __init__(self, n):
        self.Fz_list = [self.get_fz(i) for i in range(n)]


class AutoFp(Fuhrpark):
    def get_fz(self, n):
        return Auto(n)

if __name__ == '__main__':

    my_autos = AutoFp(2)

    print(type(my_autos.Fz_list[0]))










