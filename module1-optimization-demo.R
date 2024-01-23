# the function
fx <- function(x) {
    return(x**3 + x**2 - x)
}
# the first derivative
fod <- function(x) {
    return (3*x**2 + 2*x - 1)
}
# the second derivative
sod <- function(x) {
    return (6*x + 2)
}

# Gradient descent iterations
grad_descent <- function(init, lr=0.1, max_iter=20){
    xold = init
    for(i in 1:max_iter) {
        xnew = xold - lr * fod(xold)
        xold = xnew
        print(paste0("At iteration ", i, ", x=", xnew))
    }
}

grad_descent(1, lr=0.1)

# Gradient descent iterations with the learning rate as in the Barzilai-Borwein method
grad_descent_bb <- function(init, lr=0.1, max_iter=20){
    xold = init
    for(i in 1:max_iter) {
        xnew = xold - lr * fod(xold)
        lr = abs((xnew-xold) * (fod(xnew)-fod(xold))) / (fod(xnew)-fod(xold))**2
        xold = xnew
        print(paste0("At iteration ", i, ", x=", xnew))
    }
}

grad_descent_bb(1, lr=0.1)

# Newton's method
newton <- function(init, max_iter=5) {
    xold = init
    for(i in 1:max_iter) {
        xnew = xold - fod(xold)/sod(xold)
        xold = xnew
        print(paste0("At iteration ", i, ", x=", xnew))
    }
}

newton(1)
