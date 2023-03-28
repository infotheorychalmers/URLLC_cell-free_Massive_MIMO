function z = pagetimes(x,y)
% Pagetimes perform pagewise matrix product for two 3D arrays
% Input: x,y are 3D arrays (being vectors or matrix in pages)
% Output: z is the pagewise product of x and y
% 
% Zhang Tongyi Nov.22,2020

[mx, nx, p]=size(x);
[ny, py, q]=size(y);
if p~=q
    if p==1
        x = repmat(x, 1, 1, q);
    end
    if q==1
        y = repmat(y, 1, 1, p);
    end
    if p~=1 && q~=1
        error('Pagetimes Error: pages of "x" and "y" are not matched!')
    end
end
if ny==nx
    x1 = reshape(x, mx, nx, 1, size(x,3));
    y1 = reshape(y, ny, 1, py, size(x,3));
    y2 = permute(y1,[2,1,3,4]);
    z = reshape(sum(bsxfun(@times,x1,y2),2),mx,py,size(x,3));
elseif ny~=nx
    if mx==1 && nx==1
        x1 = repmat(x, ny, 1, 1);
        z = x1.*y;
    elseif ny==1 && py==1
        y1 = repmat(y, 1, nx, 1);
        z = x.*y1;
    else
        error('Pagetimes Error: dimensions of "x" and "y" are not matched for matrix times!')
    end
end
end
