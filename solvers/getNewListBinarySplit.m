function newList=getNewListBinarySplit(prevList,d)
%given a list of cubes, by subdividing each one of the cubes into two
if isempty(prevList)
    newList=prevList;
    return
end

prevListSize=length(prevList);
newList=cell(2*prevListSize,1);
listInd=0;
for ii=1:prevListSize
    cube=prevList{ii};
    new_h=cube.h;
    x0=cube.x;
    [hSplit,indSplit]=max(new_h);
    new_h(indSplit)=hSplit/2;
    B=zeros(2,d);
    B(1,indSplit)=1;
    B(2,indSplit)=-1;
    for jj=1:2
        listInd=listInd+1;
        new_cube.h=new_h;
        new_cube.x=x0+(hSplit/2)*B(jj,:)';
        newList{listInd}=new_cube;
    end
end



end



