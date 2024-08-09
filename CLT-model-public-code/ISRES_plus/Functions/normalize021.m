function y = normalize021(x)

    minx = min(x,[],'all');
    maxx = max(x,[],'all');
    
    y    = (x - minx)/(maxx - minx);
end

