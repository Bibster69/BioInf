import itertools

if __name__ == '__main__':
    list = []
    factorial = 1
    for i in range(1, 6):
        factorial = factorial * i
        list.append(i)

    print(factorial)
    for perm in itertools.permutations(list):
        strr = ''
        for p in perm:
            tmp = str(p)
            strr += tmp + ' '
        print(strr)


