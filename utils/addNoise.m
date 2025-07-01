function [N_noisy,Variance] = addNoise(Data,Noise_Ratio,Type)
signal_power = rms(Data,1);
if isequal(Type,"normal")
    sigma = Noise_Ratio.*signal_power;
    noise = normrnd(0,ones(size(Data(:,1)))*sigma,size(Data));
    N_noisy = Data + noise;
    Variance = mean(sigma);
elseif isequal(Type,"lognormal")
    sigma = Noise_Ratio;
    noise = normrnd(0,sigma,size(Data));
    N_noisy = Data .* exp(noise);
    Variance = (exp(Noise_Ratio^2)-1)*exp(Noise_Ratio^2)+(exp(Noise_Ratio^2/2)-1)^2;
else
    error("Unknown Noise Type")
end
end