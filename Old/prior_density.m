function p = prior_density(params, type, prior_params)
    switch lower(type)
        case 'gaussian'
            mu = prior_params.mean;
            sigma = prior_params.std;
            p = prod(normpdf(params, mu, sigma));

        case 'exponential'
            lambda = 1 ./ prior_params.mean;
            p = prod(exppdf(params, 1 ./ lambda));

        case 'lognormal'
            mu = log(prior_params.mean.^2 ./ sqrt(prior_params.std.^2 + prior_params.mean.^2));
            sigma = sqrt(log(1 + (prior_params.std.^2) ./ (prior_params.mean.^2)));
            p = prod(lognpdf(params, mu, sigma));

        case 'laplace'
            mu = prior_params.mean;
            b = prior_params.std / sqrt(2);
            p = prod(1/(2*b) * exp(-abs(params - mu)/b));

        case 'uniform'
            lb = prior_params.lb;
            ub = prior_params.ub;
            in_bounds = all(params >= lb & params <= ub);
            p = in_bounds * prod(1 ./ (ub - lb));

        otherwise
            error('Unknown prior type.');
    end
end
