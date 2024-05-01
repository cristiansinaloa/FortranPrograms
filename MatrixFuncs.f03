Module MatrixFuncs

CONTAINS

    Subroutine PrintMat(A,n,m)
      implicit none
      real(kind=8),intent(in),dimension(:,:) :: A
      integer,intent(in),optional :: n,m
      integer :: i,j
    
      If (present(n) .and. present(m)) then
        Do i=1,n
          write(*,'("|")',advance='no')
          Do j=1,m
            write(*,'(f8.3,t3)',advance='no') A(i,j)
          End Do
          write(*,'(" |")')
        End Do
      Else
        Do i=1,size(A,1)
          write(*,'("|")',advance='no')
          Do j=1,size(A,2)
            write(*,'(f8.3,t3)',advance='no') A(i,j)
          End Do
          write(*,'(" |")')
        End Do
      End If
    End Subroutine PrintMat

    Subroutine PrintVec(v,n)
      implicit none
      real(kind=8),intent(in),dimension(:) :: v
      integer,intent(in),optional :: n
      integer :: i
   
      If (Present(n)) then
        Do i=1,n
            write(*,'("|",f8.3,"| ")') v(i)
        End Do
      Else
        Do i=1,size(v)
            write(*,'("|",f8.3,"| ")') v(i)
        End Do
      End If
    End Subroutine PrintVec

    Subroutine slicef(Minor,A,n,row,col)
      implicit none
      integer,intent(in) :: n,row,col
      real(kind=8),dimension(n,n),intent(in) :: A
      real(8),dimension(n-1,n-1),intent(out) :: Minor
      logical,dimension(n,n) :: mask
  
      mask = .True.
      mask(row,:) = .False.
      mask(:,col) = .False.
  
      Minor = reshape(pack(A,mask),[n-1,n-1])
    End Subroutine slicef

    Recursive Function det(A,n) result(ans)
      implicit none
      real(kind=8),intent(in),dimension(n,n) :: A
      integer,intent(in) :: n  
      real(kind=8) :: ans
      real(kind=8),dimension(n-1,n-1) :: sl
      integer :: i
  
      ans = 0
  
      If (n == 2) then
        ans = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        return
      else if (n == 1) then
        ans = A(1,1)
        return
      else
        Do i = 1,n
          call slicef(sl,A,n,1,i)
          ans = ans + ((-1.0)**(1+i) * A(1,i)*det(sl,n-1))
        End Do
        return
      End if
    End Function det

    Subroutine Cram(x,A,b,n)
      implicit none
      integer,intent(in) :: n
      real(8),dimension(n,n,n+1),intent(inout) :: A
      real(8),dimension(n),intent(in) :: b
      real(8),dimension(n,1),intent(out) :: x
      real(8) :: det_A
      integer :: i
  
      det_A = det(A(:,:,1),n)
  
      Do i=1,n
        A(:,:,i+1) = A(:,:,i)
        A(:,i,i+1) = b
        x(i,1) = det(A(:,:,i+1),n) / det_A
      End Do
    End Subroutine Cram

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

End Module MatrixFuncs
