x = linspace(-0.075, 0.075, 1500);


y = zeros(size(x));

figure(1)
clf

for i=1:length(x)
    X = x(i);
    if X < -0.0275 || X > 0.0275
          y(i)=(1E9 .* X .^ 10) + 0.039;
      end
      if X >= -0.0275 && X <= 0.0275
          y(i) = ((-1 * ((0.039^5) / (X.^2 + 0.039^(4))))) + 0.039;
      end  
end








plot(x, y)
axis equal
z = zeros(size(x));

data = table(x', y', z');

writetable(data, "data.txt")