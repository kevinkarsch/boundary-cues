function names = MIT_TRAIN_EXPAND(tags)

names = MIT_TRAIN;
N = {};
for ti = 1:length(tags)
  N{ti} = cellfun(@(x) [tags{ti}, x], names, 'UniformOutput', false);
end
names = cat(2,N{:});
