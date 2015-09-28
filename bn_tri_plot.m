figure(1);hold on;
tmp;
ntri = size(tris, 1)/4;
% indx = 1;
% for i = 1:ntri
%     for j = 1:3
%         PT = [tris(indx, :);tris(indx+1, :)];
%         line(PT(:,1), PT(:,2), PT(:,3));
%         hold on;
%         indx = indx + 1;
%     end
%     indx = indx + 1;
% end

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
hold on
plot3(x(:,1), x(:,2), x(:, 3), 'k.');
