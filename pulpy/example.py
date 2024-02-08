def add_one(number):
    return number + 1

def fib(number):    
    if number <= 1:
        return number
    else:
        return fib(number-1) + fib(number-2)
