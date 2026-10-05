import itertools
import unittest
import numpy as np
from scipy.stats import binom

from weak_bias_precision import mixture_chi_square, local_length_lower_bound


class WeakBiasPrecisionTests(unittest.TestCase):
    def test_exact_vector_likelihood_moment(self):
        for size in [2, 4, 8]:
            draws = np.array(list(itertools.product([0, 1], repeat=size)))
            for shift in [.01, .04]:
                w = .8
                probabilities = [.5+shift, .5-w/(1-w)*shift]
                mixture = sum(weight*np.prod(np.where(draws, p, 1-p), axis=1)
                              for weight,p in zip([w, 1-w], probabilities))
                null = np.full(len(draws), .5**size)
                exact = np.sum((mixture-null)**2/null)
                self.assertAlmostEqual(mixture_chi_square(size, shift), exact, places=13)

    def test_count_conditioning_and_target_joint_bound(self):
        n, count, w, delta = 2, 4, .8, .04
        vectors = np.array(list(itertools.product(range(n+1), repeat=count)))
        null = np.prod(binom.pmf(vectors, n, .5), axis=1)
        alternative = np.zeros(len(vectors))
        q = binom.cdf(2, count, w)
        for flags in itertools.product([0, 1], repeat=count):
            flags = np.array(flags)
            if flags.sum() < 3:
                continue
            weight = np.prod(np.where(flags, w, 1-w))/(1-q)
            probabilities = np.where(flags, .5+delta, .5-w/(1-w)*delta)
            alternative += weight*np.prod(binom.pmf(vectors, n, probabilities), axis=1)
        np.testing.assert_allclose(alternative.sum(), 1., atol=1e-14)
        target_null, target_alt = [binom.pmf(np.arange(n+1), n, p) for p in [.5,.5+delta]]
        total_variation = np.abs(np.outer(alternative,target_alt)-np.outer(null,target_null)).sum()/2
        chi = mixture_chi_square(n, delta)
        upper = .5*np.sqrt((1+chi)**count*(1+4*delta**2)**n-1)+q
        self.assertLessEqual(total_variation, upper)

    def test_iid_randomized_full_observation_likelihood(self):
        # Four observed outcomes (A,Y); X is ancillary and cancels from ratios.
        for size in [2,4]:
            draws=np.array(list(itertools.product(range(4),repeat=size)))
            shift,w=.03,.8
            def mass(theta):
                probabilities=np.array([.3,.2,.5*(1-theta),.5*theta])
                return np.prod(probabilities[draws],axis=1)
            null=mass(.5)
            alternative=w*mass(.5+shift)+(1-w)*mass(.5-w/(1-w)*shift)
            exact=np.sum((alternative-null)**2/null)
            self.assertAlmostEqual(mixture_chi_square(size,shift,w,2.),exact,places=13)
        value=local_length_lower_bound(1000,1000,4096,experiment='iid_randomized')
        self.assertGreater(value['scaled_lower_bound'],.08)

    def test_fourth_root_scaling_and_large_target_cap(self):
        values = [local_length_lower_bound(500,500,k) for k in [256,4096,65536]]
        for value in values:
            self.assertGreater(value['scaled_lower_bound'], .07)
        self.assertAlmostEqual(values[0]['shift']/values[1]['shift'], 2.)
        large_target = local_length_lower_bound(500,10000000,256)
        self.assertAlmostEqual(large_target['shift'], .1/np.sqrt(10000000))

    def test_grouped_intercept_statistic_matches_patient_fit(self):
        from continuous_balance import linear_balance_design, evaluate_linear_balance
        from run_weak_bias_precision import intercept_summary
        design=linear_balance_design(np.ones((8,1)),np.ones(1))
        for successes in range(9):
            outcomes=np.r_[np.ones(successes),np.zeros(8-successes)]
            raw=evaluate_linear_balance(design,outcomes)
            grouped=intercept_summary(successes,8,design)
            for key in raw:
                np.testing.assert_allclose(raw[key],grouped[key],rtol=0,atol=1e-14)


if __name__ == '__main__':
    unittest.main()
