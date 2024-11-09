function result = MA_filter(signal, order)
    result = zeros(size(signal));
    
    for n = order:length(signal)
        result(n) = (1/order) * sum(signal(n-order+1:n));
    end
    
    result(1: order-1) = result(order);
end

