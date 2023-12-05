module support_subs
contains

subroutine partition(A, marker)
    implicit none
!   integer, intent(in) :: n
!   real, dimension(:), intent(inout) :: A
    real, intent(inout) :: A(:)
    integer, intent(out) :: marker

    integer :: l, h
    real    :: pk  ! pivotkey

    l = 1
    h = size(A)  ! n
    pk = A(l)

    do while (l < h)
        do while (l<h .and. A(h)>=pk)
            h = h - 1
        end do
        A(l) = A(h)

        do while (l<h .and. A(l)<=pk)
            l = l + 1
        end do
        A(h) = A(l)
    end do

    marker = l
    A(marker) = pk
end subroutine partition


recursive subroutine qsort(A)
    implicit none
!   integer, intent(in) :: n
    real, intent(inout) :: A(:)

    integer :: marker

    if (size(A) <= 1) then  ! n
        return  ! you need to have stop condition for the recursion.
    else
        call partition(A, marker)
        call qsort(A(:marker-1))  ! Ligang: what is marker is already 1 or size(A)? segment fault? 
        call qsort(A(marker+1:))  !   or size(A) == 1?
    end if
end subroutine qsort

end module support_subs
