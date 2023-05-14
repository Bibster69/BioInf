

def fibo(n, lifespan):

    gen = 1
    parents = 1
    kids = 0
    while gen <= n - 1:
        sum = parents + kids
        kids = parents
        parents = sum
        gen += 1

    print(parents)


if __name__ == '__main__':
    print(fib(6, 3))
