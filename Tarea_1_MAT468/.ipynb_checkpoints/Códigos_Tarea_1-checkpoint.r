# Datos
A = matrix(c(1,2,1,1,3,-1,1,1,-1),3,3)
b = matrix(c(4,9,-1),3,1)

# Descomposición LU
LU = function( A ){
    if ( dim(A)[1]!=dim(A)[2] ){
        stop( 'La matriz no es cuadrada' )
    }
    if ( !is.numeric( A ) ){
        stop( 'Argumentos Erroneos' )
    }
    n = dim(A)[1]
    L = matrix(0,n,n)
    U = diag(n)
    for ( i in 1:n ) {
    aux1 = i - 1 
    for (j in 1:n){
        L[j,i] = A[j,i]
        if ( aux1 > 0 ) {
                for ( k in 1:aux1 ) {
                    L[j,i] = L[j,i] - L[j,k] * U[k,i]
                }
            }  
    }
    aux2 = i+1 
    if ( aux2 <= n ) {
            for ( j in aux2:n ) {
                U[i,j] = A[i,j] 
                if ( aux1 > 0 ) {
                    for ( k in 1:aux1 ) {
                        U[i,j] = U[i,j] - L[i,k] * U[k,j]
                    }
                }
                U[i,j] = U[i,j] / L[i,i]
            }    
        }
    }
    return(list(L=L,U=U))
}

LU(A)

#Resolver triangular inferior
rti = function (L,b){
    n = dim(L)[1]
    z = c()
    for (i in 1:n){
        sol = b[i]
        if ( length(z) > 0 ) {
                for ( k in 1:length(z) ) {
                    sol = sol - L[i,k]*z[k]
                }
            }
        sol = sol/L[i,i]
        z = c(z,sol)
    }
    return(matrix(z,n,1))
}

#Resolver triangular superior
rts = function (U,b){
    n = dim(U)[1]
    z = c()
    for (i in n:1){
        sol = b[i]
        if ( length(z) > 0 ) {
                for ( k in 1:length(z) ) {
                    sol = sol - U[i,n-k+1]*z[k]
                }
            }
        sol = sol/U[i,i]
        z = c(sol,z)
    }
    return(matrix(z,n,1))
}

#Resolver un sistema Ax = b por LU
resolver = function (A,b){
    L =LU(A)$L
    U =LU(A)$U
    z=rti(L,b)
    x = rts(U,z)
    return(x)
    }

#comparar con la función del R
resolver(A,b)
solve(A,b)

# Refinamiento iterativo
refinamiento_iterativo = function (kmax,tol,A,b,x_0){
    k=0
    while (k<kmax){
        r = b - A%*%x_0
        delta = resolver(A,r)
        x_0 = x_0  + delta
        if (max(delta)<= tol*max(x_0)){
            return(x_0)
        }
        k=k+1
    }
        stop(paste('El algoritmo no converge despues de',kmax,'iteraciones'))
}

refinamiento_iterativo(10,0.001,A,b,rnorm(3))

# datos
B = matrix(c(30,16,46,16,10,26,46,26,72),3,3)

# Función sweep
sweep = function(A, k) {
    B = matrix(0,dim(A)[1],dim(A)[1])
    B[k,k] = 1/A[k,k]
    for (i in 1:n){
        if (i != k){
            B[i,k]=-A[i,k]/A[k,k]
        }
    }
    for (j in 1:n){
        if (j != k){
            B[k,j]=A[k,j]/A[k,k]
        }
    }
    for (i in 1:n){
        for (j in 1:n){
            if (i != k & j != k){
                B[i,j]= A[i,j]-(A[i,k]*A[k,j])/A[k,k]
            }
        }
    }
    return(B)
}

#Invertir mediante sweep
inv_sweep = function (A) {
    n = dim(A)[1]
    B = A
    for (i in 1:n){
        B = sweep(B,i)
    }
    return(B)
}

# Comparar con una función del R
inv_sweep(A)
solve(A)

# Función Cholesky
choleski = function(A){
    n = dim(A)[1]
    T = matrix(0,n,n)
    T[1,1] = sqrt(A[1,1])
    for (j in 2:n){
        T[1,j] = A[1,j]/T[1,1]
    }
    for (i in 2:n){
        res = A[i,i]
        for (k in 1:(i-1)){
            res = res - (T[k,i])^2
        }
        T[i,i] = sqrt(res)
        res=0
        if ((i+1)<=n){
        for (j in (i+1):n){
            res = A[i,j]
            for (k in 1:(i-1)){
                res = res - T[k,i]*T[k,j]
            }
            T[i,j] = res/T[i,i]
            }
        }
    }
    return(T)
}

choleski(B)

# Comentarios
# multiplicar matrices
multiplicar_matrices = function (A,B){
    if (dim(A)[2]!=dim(B)[1]){
        stop('No coinciden las dimensiones')
    }
    m = dim(A)[1]
    n = dim(A)[2]
    p = dim(B)[2]
    C = matrix(,m,p)
    for (i in 1:m){
        for (j in 1:p){
            C[i,j]=0
            for (k in 1:n){
                C[i,j] = C[i,j]+ A[i,k]*B[k,j]
            }
        }
    }
    return(C)
}

# Descomposició LU versión 2
LU_ver2 = function( A ){
    if ( dim(A)[1]!=dim(A)[2] ){
        stop( 'La matriz no es cuadrada' )
    }
    if ( !is.numeric( A ) ){
        stop( 'Argumentos Erroneos' )
    }
    n = dim(A)[1]
    L = matrix(0,n,n)
    U = diag(n)
    diag(L) = rep(1,n)
    for (i in 1:n) {
        aux1 = i - 1
        for (j in 1:n) {
            U[i,j] = A[i,j]
            if ( aux1 > 0 ) {
                for ( k in 1:aux1 ) {
                    U[i,j] = U[i,j] - L[i,k] * U[k,j]
                }
            }
        }
        aux2 = i + 1
        if ( aux2 <= n ) {
            for ( j in aux2:n ) {
                L[j,i] = A[j,i]
                if ( aux1 > 0 ) {
                    for ( k in 1:aux1 ) {
                        L[j,i] = L[j,i] - L[j,k] * U[k,i]
                    }
                }
                L[j,i] = L[j,i] / U[i,i]
            }    
        }
    }
    return(list(L=L,U=U))
}
