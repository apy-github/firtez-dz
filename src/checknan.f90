!
! This function checks whether the argument is a valid floating point number
Logical Function CheckNaN(Number)
!
Implicit None
Real(8) :: Number
!
!print*,Number,'I am NAN checker'
CheckNaN=.False.
If (Number .gt. Huge(Number)) then 
!   Print *,'Larger than Huge'
   CheckNan=.True.
End if
If (Number .lt. -Huge(Number)) then
!   Print *,'Smaller than -Huge'
   CheckNan=.True.
End if
If (Number .ne. Number) then
!   Print *,'Not equal to itself'
   CheckNan=.True.
End if
If ( (Number .gt. 0.0) .EQV. (Number .le. 0.0) ) then
!   Print *,'Does not pass the EQV test'
   CheckNan=.True.
End if
If (Abs(Number) .lt. 0.) then
!   Print *,'Negative absolute value'
   CheckNan=.True.
End if
!
if (Number.eq.(Number+1d0)) then
!   Print *,'Value eq to value', Number, (Number+1d0)
   CheckNan=.True.
End if
Return
End Function CheckNaN
!
