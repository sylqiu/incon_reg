function pind = findclstpt(V,p)
pind = zeros(size(p,1),1);
for i = 1:size(p,1)
    x = p(i,:);
    [~,ind] = min(sum((V - repmat(x,size(V,1),1)).^2,2));
    pind(i) = ind; 
end
end