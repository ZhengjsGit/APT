function idx=GAPS_APT_TensorPointer(Index,Dim,Order)
    if length(Index)~=Order
        error('The length of "Index" does NOT equal to "Order"');
        idx=0;
    else
        if 0==Order
            idx=1;
        else
            idx=Index(Order)-1;
            for i=Order-1:-1:1
                idx=idx*Dim+(Index(i)-1);
            end
            idx=idx+1;
        end
    end
end