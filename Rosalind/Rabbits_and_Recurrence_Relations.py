

def fibo(n, k):

    gen = 1
    parents = 1
    kids = 0
    while gen <= n - 1:
        sum = parents + kids
        kids = parents * k
        parents = sum
        gen += 1


    print(parents)

def fibo_rec(n, k):
    if n == 1:
        return 1
    if n == 2:
        return k

    if n == 3:
        return fibo_rec(n-1, k) + fibo_rec(n-2, k)

    if n == 4:
        return fibo_rec(n - 1, k) + fibo_rec(n - 2, k)

    return fibo_rec(n-1, k) + (fibo_rec(n - 2, k) * k) # We can see that the number of rabbits two generations back is equal to the number of adults one generation back.



if __name__ == '__main__':
    print(fibo(33, 3))
    print(fibo_rec(5, 3))