def multiply(x, y):
    return (x * y)


def get5multiplier():
    return lambda x: multiply(5, x)


const = get5multiplier()

print(const(4))