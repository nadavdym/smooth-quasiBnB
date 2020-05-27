function boundaryList=getBoundaryList(preBoundaryList,c0)
%given list of cubes, each of d-dimension, return all d-1 faces which are
%common with c0
if isempty(preBoundaryList)
    boundaryList=preBoundaryList;
    return
end
boundaryList={};
for ii=1:length(preBoundaryList)
    c=preBoundaryList{ii};
    d=length(c.h);
    non_zero_coords=find(c.h);
    nnz=length(non_zero_coords);
    for ii=1:nnz
        ind=non_zero_coords(ii);
        if c0.x(ii)+c0.h(ii)==c.x(ind)+c.h(ind)
            new_c=c;
            new_c.h(ind)=0;
            new_c.x(ind)=c.x(ind)+c.h(ind);
            boundaryList={boundaryList{:}, new_c};
        end
        if c0.x(ii)-c0.h(ii)==c.x(ii)-c.h(ii)
            new_c=c;
            new_c.h(ind)=0;
            new_c.x(ind)=c.x(ind)-c.h(ind);
            boundaryList={boundaryList{:}, new_c};
        end
    end

end