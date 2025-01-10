import sys
import math

def fibonacci(a, b, n):
    ''' The goal of this code is to calculate the n-th term of a sequence without generating every term up to n. 
        This approach avoids recursion or loops by directly computing the result with an explicit formula.
    '''
        
    # eigenvalues: roots of the characteristic equation for the Fib sequence (represent the growth rates)
    phi = (1 + math.sqrt(5)) / 2 # golden ratio (φ) 
    psi = (1 - math.sqrt(5)) / 2 # conjugate (ψ)
    
    # we solve for the coefficients c1 and c2 using the initial values a and b. 
    # these constants adjust the formula to match the specific starting values given by the user
    c1 = (b - psi * a) / (phi - psi) 
    c2 = (phi * a - b) / (phi - psi)
    
    # each term in the sequence can be expressed as a combination of powers of the eigenvalues
    # by using linear combinations of these terms, we can represent any term as the sum of two terms:
    Fn = c1 * (phi ** (n-1)) + c2 * (psi ** (n-1)) # adjust the formula to the n-th term
    
    # return the integer value of the nth term
    return int(round(Fn))


if __name__ == "__main__":
    a = int(sys.argv[1])
    b = int(sys.argv[2])
    n = int(sys.argv[3])
    result = fibonacci(a, b, n) # compute the nth term
    print(result)


'''in a standard Fibonacci sequence, each term is the sum of the two preceding terms:
''' # F(n) = F(n-1) + F(n-2)
'''
matrix representation:
    we can represent this recurrence relation in matrix form as follows:
''' # | F(n)   |   =   |1  1|   *   |F(n-1)|
    # | F(n-1) |       |1  0|       |F(n-2)|
    
    # first row: F(n)=1⋅F(n−1)+1⋅F(n−2)=F(n−1)+F(n−2) => F(n) = F(n-1) + F(n-2)
    # second row: (n−1)=1⋅F(n−1)+0⋅F(n−2)=F(n−1) => F(n-1) = F(n-1)
'''    
    each time you apply the matrix multiplication, you calculate the next two terms in the Fibonacci sequence.

repeated multiplication:
    using this matrix form, each step in the Fibonacci sequence can be computed
    by repeatedly multiplying this transformation matrix.
    We can express the n-th term as:
''' # | F(n)   | = |1  1|^(n-1) * | F(2) |
    # | F(n-1) |   |1  0|         | F(1) |
'''
efficient calculation:
    this matrix form allows us to compute F(n) directly by calculating
    the power of the matrix A = |1 1; 1 0| raised to (n-1).
    this avoids the need for recursive calculations and gives us a direct way to find any term F(n).

eigenvalues:
    constants (φ and ψ) are the roots of the characteristic equation of this matrix, derived from: det(A - lambda * I) = 0
    by raising the eigenvalues to powers, we capture the growth patterns of the sequence
    These roots capture the seq's growth behavior and allow us to express the Fibonacci seq in a closed-form solution.
'''