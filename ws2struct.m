function p = ws2struct
% pack the whole workspace into a structure
% O.Codol 22nd Jan. 2020
% codol.olivier@gmail.com
%---------------------------------------------

varnames = evalin('caller', 'who');
for k = 1:size(varnames, 1)
    thisval = evalin('caller', varnames{k});
    p.(varnames{k}) = thisval;
end

end

