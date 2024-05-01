Module NewtonRaphson2D_mod
    implicit none
    
CONTAINS

    Subroutine NewtonRaphson2D(f, x0, x_min,tol)
        ! This subroutine finds the minimum of a function using the Newton Raphson method.
        ! Inputs:
        !   f: Function to minimize (must be a function of one variable)
        !   x0: Initial guess for the minimum
        ! Outputs:
        !   x_min: The minimum of the function
        ! Optional Inputs: 
        !   tol: Tolerance for convergence
        
        implicit none
        Interface
            Function f(x)
                real :: f
                real, intent(in),dimension(:) :: x
            End Function f
        End Interface
        real, intent(in),dimension(:) :: x0
        real, optional :: tol
        real,intent(out),dimension(size(x0)) :: x_min
        real,dimension(size(x0)) :: xk,g
        real,dimension(size(x0),size(x0)) :: H,Hinv
        real :: def_tol
        integer :: k,maxiters

        If (present(tol)) Then
            def_tol=tol
        Else
            def_tol = 1e-8 ! max possible dx since round off error kicks in
                           ! for finite differences of the gradient and hessian
        End If

        xk = x0
        write(*,*) "====== Please Wait Locating Minimum... ======"

        maxiters = 500
        k = 0
        Do While(k < maxiters)
            ! Calculate the gradient of the function f at the pt xk
            g = gradient(f,xk,def_tol)

            ! Check gradient convergence before Hessian Calc
            If (abs(norm2(g)) < def_tol) Exit

            ! Alternative gradient convergence criteria using Infinity norm
            !If (maxval(abs(g)) < def_tol) Exit

            H = Hessian(f,xk,def_tol)
            Hinv = inv(H)
            
            ! Update x using Newton's method
            xk = xk - matmul(Hinv,g)

            ! Print statement for debugging and seeing the iterates
            !write(*,*) xk

            k = k + 1
        End Do
        write(*,*) "====== N-Dim Newton Raphson Method for Minimization ======"
        write(*,'(X,A,ES12.2)') "tolerance =", tol
        write(*,*) "|grad(f)| =", abs(norm2(g))
        write(*,*) "Iterations =", k
        write(*,*) "x-val =", xk
        write(*,*) "f-val = ", f(xk)
        
        x_min = xk
        
    End Subroutine NewtonRaphson2D
    
    Function gradient(f,x0,tol) result(g)
        ! Function to calculate the gradient of a function using centered finite differences
        ! Inputs:
        !   f(function): Function to minimize
        !   x0(Vector): Evaluate the gradient at that point
        ! Outputs:
        !   g(Vector): The gradient of f evaluated at x0
        implicit none
        Interface
            Function f(x)
                real :: f
                real, dimension(:),intent(in) :: x
            End Function f
        End Interface
        real,intent(in),dimension(:) :: x0
        real,intent(in),optional :: tol
        real,dimension(size(x0)) :: g
        real,dimension(size(x0)) :: ei
        real :: dx
        integer :: i,n

        If (present(tol)) Then
            dx=tol
        Else
            dx = 1e-8 ! max possible dx since round off error kicks 
                      ! in for finite difference
        End If

        n = size(x0)
        g = 0

        Do i=1,n
            ei = 0
            ei(i) = 1
            g(i) = (f(x0 + dx*ei) - f(x0 - dx*ei)) / (2.0 * dx)
        End Do
        
    End Function gradient

    Function Hessian(f,x0,tol) result(H)
        ! Function to calculate the Hessian of a vector function using centered finite differences
        ! Inputs:
        !   f(function): Function to compute the Hessian of 
        !   x0(Vector): Evaluate the Hessian at that point
        ! Outputs:
        !   H(Matrix): nxn Hessian Matrix 
        implicit none
        Interface
            Function f(x)
                real :: f
                real, intent(in),dimension(:) :: x
            End Function f
        End Interface
        real,intent(in),dimension(:) :: x0
        real,intent(in),optional :: tol
        real,dimension(size(x0),size(x0)) :: H
        real,dimension(size(x0)) :: ei,ej
        real :: dx
        real :: f1,f2,f3,f4
        integer :: i,j,n

        If (present(tol)) Then
            dx = tol
        Else
            dx = 1e-8 ! max possible dx since round off error kicks 
                      ! in for finite difference
        End If

        n = size(x0)
        H = 0

        Do i=1,n
            Do j=1,n
                ei = 0
                ei(i) = 1
                ej = 0
                ej(j) = 1
                f1 = f(x0 + dx*ei + dx*ej)
                f2 = f(x0 + dx*ei - dx*ej)
                f3 = f(x0 - dx*ei + dx*ej)
                f4 = f(x0 - dx*ei - dx*ej)
                H(i,j) = (f1-f2-f3+f4) / (4*dx*dx)
            End Do
        End Do
        
    End Function Hessian

    Function inv(A) result(Ainv)
        ! Function to calculate the inverse of matrix using gaussian elimination
        ! Inputs:
        !   A(Array): NxN Matrix
        ! Outputs:
        !   Ainv(Array): The inverse of the input matrix
        implicit none
        real, dimension(:,:) :: A
        real, dimension(size(A,1), size(A,2)) :: Ainv
        real :: coeff
        integer :: i,j,n
        logical :: singularity
        
        n = size(A,1)

        ! Set Ainv equal to the identity matrix 1st
        Ainv = 0
        Do i=1,n
            Ainv(i,i) = 1
        End Do
        
        Do i=1,n
            ! Check for singularity
            If (A(i,i) == 0.0) Then
                singularity = .True.
                write(*,*) "Matrix is singular cannot invert"
                Exit
            End If
            
            ! Scale row i so A(i,i) becomes 1
            coeff = A(i,i)
            A(i,:) = A(i,:) / coeff
            Ainv(i,:) = Ainv(i,:) / coeff
            
            ! Subtract multiples of row i from other rows to make A(:,i) become 0
            Do j=1,n
                If (i /= j) Then
                    coeff = A(j,i) / A(i,i)
                    A(j,:) = A(j,:) - coeff * A(i,:)
                    Ainv(j,:) = Ainv(j,:) - coeff * Ainv(i,:)
                EndIf
            End Do

        End Do
        
    End Function inv

End Module NewtonRaphson2D_mod

Program NewtonRaphson2DTest
    ! compile as: gfortran -O3 -fdefault-real-8 NR2D.f03 -o NR2D.exe
    ! this gives the best accuracy and performance
    use NewtonRaphson2D_mod
    implicit none
    real,dimension(2) :: x_min,x_guess
    real :: tic, toc,tol
    ! real,dimension(2,2) :: B, Binv 
    ! real,dimension(3,3) :: Hcheck
    ! real,dimension(3) :: xtest,gcheck
    ! integer :: i,j

    x_guess = [0.0,0.0]
    tol = 1e-8
    ! Find the minimum of the function using Newton's method
    call cpu_time(tic)
    call NewtonRaphson2D(f,x_guess,x_min,tol)
    call cpu_time(toc)

    write(*,*) "Minimum of the function at x =", x_min

    write(*,'(X,A,ES10.2,A)') "time = ", toc-tic, " sec"

    ! xtest = [1.0,1.0,1.0]
    ! gcheck = gradient(f,xtest)

    ! Do i=1,3
    !     write(*,*) gcheck(i)
    ! End Do

    ! xtest = [1.0,1.0,1.0]
    ! Hcheck = Hessian(f,xtest)

    ! Do i=1,3
    !     write(*,*) (Hcheck(i,j), j=1,3)
    ! End Do


    ! This is code to test the inv function the inverse is 
    ! Ainv = [[2,3],[2,2]]
    ! B = transpose(reshape([-1.,1.5,1.,-1.],[2,2]))
    ! Binv = inv(B)

    ! Do i=1,2
    !     write(*,*) (Binv(i,j), j=1,2)
    ! End Do
    
    
    CONTAINS
    ! Define the function to minimize
    Function f(x)
        real :: f
        real,intent(in),dimension(:) :: x

        ! Rosenbrock function for testing x_min @ [1,1]
        f = (1.0 - x(1))**2 + 100.0*(x(2) - x(1)**2)**2 

        ! Multimodal function x_min @ [3.0,0.5]
        !f = (1.5 - x(1) + x(1)*x(2))**2 + (2.25 - x(1) + x(1)*x(2)**2)**2 + (2.625 - x(1) + x(1)*x(2)**3)**2

        ! gradient Test Func w/ vec [1,1,1] ans = [4,4,2]
        !f = x(1)**2 + x(2)**2 + 2*x(1)*x(2)*x(3)

        ! Hessian Test Func w/ vec [1,1,1] ans = [[2,2,1],[2,2,1],[1,1,0]]
        !f = x(1)**2 + x(2)**2 + x(1)*x(2) + x(3) + x(1)*x(2)*x(3)

        ! Example function: x^2 - 4x + 3; x_min @ x= 2
        !f = x(1)**2 - 4*x(1) + 3

        ! Ex func w/ multiple minima; global x_min @ [5.145735]
        !f = sin(x(1)) + sin((10./3.)*x(1)) 
    End Function f
    
End Program NewtonRaphson2DTest