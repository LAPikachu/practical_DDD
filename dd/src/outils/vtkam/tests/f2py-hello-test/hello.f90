subroutine pow2(in_x, out_x)
  implicit none
  real, intent(in)   :: in_x
  real, intent(out)     :: out_x
  out_x = in_x ** 2
end subroutine pow2
