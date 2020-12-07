warn = warning('query','last');
if ~isempty(warn)
id=warn.identifier;
    if ~isempty(id)
        warning('off',id);
    end
end