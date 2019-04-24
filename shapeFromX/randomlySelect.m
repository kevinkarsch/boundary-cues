function X_selected = randomlySelect(X, N)

if N >= size(X,1)
  X_selected = X;

else
  
  sortidx = randperm(size(X,1));
  X_selected = X(sortidx(1:N),:,:,:);  
  
end


