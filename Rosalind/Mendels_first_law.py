
def mendels_law(k, m, n):
    recessive_odds = (1/4 * m * (m - 1) + 1/2 * m * n + 1/2 * n * m + n * (n - 1))/((k + m + n) * ((k + m + n) - 1))
    print(1 - recessive_odds)


if __name__ == '__main__':
    mendels_law(17, 15, 23)







