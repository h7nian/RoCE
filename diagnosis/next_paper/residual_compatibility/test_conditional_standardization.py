"""Patient-level independent checks of the exact-balance diagnostic."""
import itertools
import unittest
import numpy as np
from conditional_standardization import (block_dispersion_constants,block_dispersion_upper,
                                         conditional_standardization_interval)
from fitted_strata import fitted_stratum_data


def patient_matrix(coefficients,counts):
    sites,cells=counts.shape
    block=np.repeat(np.arange(sites*cells),counts.ravel())
    site=block//cells
    weights=coefficients.ravel()[block]/counts.ravel()[block]
    projection=np.eye(sites)/sites-np.ones((sites,sites))/sites**2
    matrix=projection[site[:,None],site]*weights[:,None]*weights
    for group,size in enumerate(counts.ravel()):
        index=np.ix_(block==group,block==group)
        matrix[index]*=size/(size-1)
    np.fill_diagonal(matrix,0)
    return matrix,projection,block,site,weights


class ConditionalStandardizationChecks(unittest.TestCase):
    def test_heterogeneous_block_quadratic_constants(self):
        rng=np.random.RandomState(209)
        for _ in range(50):
            counts=rng.randint(2,6,(3,4))
            coefficients=rng.normal(size=(3,4))
            probabilities=rng.uniform(.1,.9,(3,4))
            matrix,projection,blocks,sites,weights=patient_matrix(coefficients,counts)
            means=np.sum(coefficients*probabilities,axis=1)
            patient_mean=probabilities.ravel()[blocks]
            np.testing.assert_allclose(matrix@patient_mean,weights*(projection@means)[sites],atol=1e-12)
            self.assertAlmostEqual(patient_mean@matrix@patient_mean,np.var(means),places=12)
            linear,quadratic,scale,_=block_dispersion_constants(coefficients[None],counts[None])
            gaussian_matrix=matrix/4
            self.assertAlmostEqual(quadratic[0],2*np.sum(gaussian_matrix**2),places=12)
            self.assertAlmostEqual(scale[0],-2*np.linalg.eigvalsh(gaussian_matrix)[0],places=12)
            vector=.5*(matrix@patient_mean)
            self.assertLessEqual(4*vector@vector,linear[0]*np.var(means)+1e-12)

    def test_exact_binary_unbiasedness_and_upper_tail_inversion(self):
        coefficients=np.array([[.3,-.6],[.3,-.6]])
        counts=np.full((2,2),2)
        probabilities=np.array([[.1,.5],[.7,.2]])
        observations=np.array(list(itertools.product([0.,1.],repeat=8))).reshape(-1,2,2,2)
        patient_probability=np.repeat(probabilities[...,None],2,axis=2)
        law=np.prod(np.where(observations==1,patient_probability,1-patient_probability),axis=(1,2,3))
        means=observations.mean(axis=3)
        estimates=np.sum(coefficients*means,axis=2)
        variances=np.sum(coefficients**2*means*(1-means),axis=2)
        dispersion=np.var(np.sum(coefficients*probabilities,axis=1))
        estimate=np.var(estimates,axis=1)-.25*variances.sum(axis=1)
        self.assertAlmostEqual(law@estimate,dispersion,places=12)
        for alpha in [.01,.05,.2]:
            upper=block_dispersion_upper(estimates,variances,np.broadcast_to(coefficients,(len(estimates),2,2)),
                np.broadcast_to(counts,(len(estimates),2,2)),alpha)
            self.assertLessEqual(law@(dispersion>upper),alpha+1e-12)

    def test_exact_calibration_identity_and_patient_budget(self):
        setting=dict(source_count=8,shared_scale=.3,pattern='weak_disjoint',seed_offset=912)
        data=fitted_stratum_data(setting,20,retain_source_counts=True)
        counts=data['source_training_counts']+data['source_held_counts']
        np.testing.assert_array_equal(counts.sum(axis=(2,3,4)),1000)
        target=data['target_training_counts']+data['target_evaluation_counts']
        masses=target.sum(axis=(2,3))/1000
        sizes=counts.sum(axis=4)
        cell_means=counts[...,1]/sizes
        prediction=data['common_prediction'][:,None,:,:]
        assisted=np.sum(masses[:,None,:,None]*prediction,axis=2)+np.sum(masses[:,None,:,None]*(cell_means-prediction),axis=2)
        result=conditional_standardization_interval(data,6)
        np.testing.assert_allclose(assisted,result['source_arm_estimates'],atol=1e-12)
        self.assertTrue(np.all((result['interval'][:,0]<=.1)&(.1<=result['interval'][:,1])))
        self.assertFalse(result['unavailable'].any())

    def test_sparse_design_fallback_and_overlap_contract(self):
        data=fitted_stratum_data(dict(source_count=4,shared_scale=.3,pattern='all_valid',seed_offset=117),2,retain_source_counts=True)
        data['source_training_counts'][0,0,0,0,:]=0
        data['source_held_counts'][0,0,0,0,:]=0
        result=conditional_standardization_interval(data,3)
        self.assertTrue(result['unavailable'][0])
        self.assertTrue(np.all(np.isfinite(result['interval'])))
        with self.assertRaises(ValueError):conditional_standardization_interval(data,2,dispersion_mode='contrast')

    def test_projection_geometry_for_finite_target_effect(self):
        data=fitted_stratum_data(dict(source_count=32,shared_scale=.3,pattern='strong',seed_offset=513),100,retain_source_counts=True)
        result=conditional_standardization_interval(data,24,point_rule='target_projection')
        counts=data['target_training_counts']+data['target_evaluation_counts']
        masses=counts.sum(axis=(2,3))/1000
        cohort_truth=masses@np.array([-.2,-.2,.4,.4])
        source=result['source_cohort_interval']
        good=(source[:,0]<=cohort_truth)&(cohort_truth<=source[:,1])
        self.assertTrue(np.any(good))
        self.assertTrue(np.all(abs(result['point'][good]-cohort_truth[good]) <=
                               abs(result['target_point'][good]-cohort_truth[good])+1e-12))
        self.assertTrue(np.all(result['point']>=result['interval'][:,0]))
        self.assertTrue(np.all(result['point']<=result['interval'][:,1]))

    def test_one_cell_reduces_to_existing_dispersion_constants(self):
        from concentrated_bounds import exponential_dispersion_constants
        coefficients=np.array([[[.3],[-.4],[.9]]])
        sizes=np.array([[[20],[30],[40]]])
        first=block_dispersion_constants(coefficients,sizes)
        second=exponential_dispersion_constants(abs(coefficients[:,:,0]),sizes[:,:,0],(1,3))
        for a,b in zip(first[:3],second):np.testing.assert_allclose(a,b,atol=1e-12)


if __name__=='__main__':unittest.main(verbosity=2)
