program integer_sum
    implicit none
    integer::number,sum
    
    sum=0
    do
        print*, "give a number(-1 to exit):"
        read*, number
        ! check for exit
        if(number==-1) then
            exit
        endif
        sum=sum+number
    enddo
    print*, "The sum is", sum
endprogram

