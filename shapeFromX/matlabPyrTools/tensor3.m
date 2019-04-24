function f = tensor3(a, b, c)

a = a(:);
b = b(:);
c = c(:);

f = bsxfun(@times, bsxfun(@times, a, b'), permute(c, [2,3,1]));
