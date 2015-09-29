figure(1);hold on;
opencascade_faces;
ntri = size(tris, 1)/4;

X = [];
Y = [];
Z = [];
for i = 1:ntri
    indx = 4*(i-1) + 1;
    X = [X, [tris(indx, 1); tris(indx+1, 1); tris(indx+2, 1)]];
    Y = [Y, [tris(indx, 2); tris(indx+1, 2); tris(indx+2, 2)]];
    Z = [Z, [tris(indx, 3); tris(indx+1, 3); tris(indx+2, 3)]];
end

C= zeros(size(X));
fill3(X,Y,Z,C);

axis equal
hold off

