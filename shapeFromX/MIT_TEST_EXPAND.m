function names = MIT_TEST_EXPAND(tags)

names = MIT_TEST;
N = {};
for ti = 1:length(tags)
  N{ti} = cellfun(@(x) [tags{ti}, x], names, 'UniformOutput', false);
end
names = cat(2,N{:});
