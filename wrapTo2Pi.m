function output = wrapTo2Pi(input)
    while input > 2*pi
        input = input - 2 * pi;
    end
    while input < 0 
        input = input + 2 * pi;
    end
    output = input;
end