x = linspace(-0.07, 0.07, 280);
a = 0.037;
n = 3.5;

y = zeros(size(x));

figure(1)
clf

for i=1:length(x)
    X = x(i);
    if (X < -0.04185 || X > 0.04185)
        y(i) = 10000000000 * X^(10) + 0.032;
        
    end
    if (X > -0.04185 && X < 0.04185)
        y(i) = (-(a^n) / (X^2 + a^(n-1))) + 0.037;
    end
end

%%parabola = (x < -0.04185 || x > 0.04185);
%%witch = (x >= -0.04185 && x < 0.04185);









plot(x, y)
axis equal
z = zeros(size(x));

data = table(x', y', z');

writetable(data, "data.txt")