function dispStruct(name, stru)

if isstruct(stru)
  for f = fieldnames(stru)'
    f = f{1};

    if isstruct(stru.(f))

      dispStruct([name, '.', f], stru.(f));

    elseif isnumeric(stru.(f))
            
      fprintf('%s.%s = \t[ ', name, f);
      fprintf('%g ', stru.(f))
      fprintf('];\n');

    else

      fprintf('%s.%s = \t{ ', name, f);
      fprintf('%g ', cell2mat(stru.(f)))
      fprintf('};\n');

    end

  end
else
  
  if iscell(stru)
    fprintf('%s = \t{ ', name);
    fprintf('%g ', cell2mat(stru));
    fprintf('}\n');
  else
    fprintf('%s = %g \n', name, stru)
  end
  
end