module MatrixFuncs

contains

  subroutine PrintMat(A,n,m)
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

  end subroutine PrintMat

  subroutine PrintVec(v,n)
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
    end subroutine PrintVec

  subroutine slicef(Minor,A,n,row,col)
    implicit none
    integer,intent(in) :: n,row,col
    real(kind=8),dimension(n,n),intent(in) :: A
    real(8),dimension(n-1,n-1),intent(out) :: Minor
    logical,dimension(n,n) :: mask

    mask = .True.
    mask(row,:) = .False.
    mask(:,col) = .False.

    Minor = reshape(pack(A,mask),[n-1,n-1])
  end subroutine slicef

  recursive function det(A,n) result(ans)
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
  end function det

  subroutine Cram(x,A,b,n)
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
  end subroutine Cram

end module MatrixFuncs
