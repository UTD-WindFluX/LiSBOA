%find external point in a 3D logical matrix
%08/09/209: created
%08/28/2019: adapted for 2D
function out=Find_bound3D(M)
    [n1,n2,n3]=size(M);
    out=zeros(size(M));
    i=repmat((1:n1)',1,n2,n3);
    j=repmat(1:n2,n1,1,n3);
    k=permute(repmat((1:n3)',1,n2,n1),[3 2 1]);
    
    face{1}=[vert(i(:,:,1)),vert(j(:,:,1)),vert(k(:,:,1))];
    face{2}=[vert(i(:,:,end)),vert(j(:,:,end)),vert(k(:,:,end))];
    face{3}=[vert(i(1,:,:)),vert(j(1,:,:)),vert(k(1,:,:))];
    face{4}=[vert(i(end,:,:)),vert(j(end,:,:)),vert(k(end,:,:))];
    face{5}=[vert(i(:,1,:)),vert(j(:,1,:)),vert(k(:,1,:))];
    face{6}=[vert(i(:,end,:)),vert(j(:,end,:)),vert(k(:,end,:))];
    
    size_k=size(unique(k));
    if size_k(1)>1
        for ID_face=1:length(face)
            B=imfill(M,face{ID_face});
            out(B>M)=1;
        end
    else
        for ID_face=3:length(face)
            B=imfill(M,face{ID_face}(:,1:2));
            out(B>M)=1;
        end
    end
    out=out==1;
end