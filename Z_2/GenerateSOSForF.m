function Fsos=GenerateSOSForF(F)
Fsos=CZprod(F{1},F{1});
for i=2:length(F)
    Fsos=CZsum(Fsos,CZprod(F{i},F{i}));
end
end