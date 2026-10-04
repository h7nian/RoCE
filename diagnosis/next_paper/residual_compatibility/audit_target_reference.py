"""Independent large-repeat check of the fitted full-target AIPW comparator."""
import json
from pathlib import Path
import sys
import numpy as np
from scipy.stats import norm


def audit(repeats=20000):
    rows=[]
    for number,heterogeneity in enumerate([0.,.3]):
        cells=np.array([(a,b) for a in [-1,1] for b in [-1,1]])
        propensity=.5+.1*cells[:,0]-.05*cells[:,1]
        outcome=np.column_stack((.55+heterogeneity*cells[:,0]/2+.08*cells[:,1],
                                 .45-heterogeneity*cells[:,0]/2+.08*cells[:,1]))
        probability=np.empty((4,2,2))
        for cell in range(4):
            for arm in range(2):
                treatment=propensity[cell] if arm==0 else 1-propensity[cell]
                probability[cell,arm]=.25*treatment*np.array([1-outcome[cell,arm],outcome[cell,arm]])
        rng=np.random.RandomState(981110+number)
        first=rng.multinomial(500,probability.ravel(),size=repeats).reshape(repeats,4,2,2)
        second=rng.multinomial(500,probability.ravel(),size=repeats).reshape(repeats,4,2,2)
        sums=[];squares=[]
        for training,evaluation in [(first,second),(second,first)]:
            arm_counts=training.sum(axis=3);cell_counts=arm_counts.sum(axis=2)
            regression=(training[:,:,:,1]+.5)/(arm_counts+1)
            treatment=(arm_counts[:,:,0]+.5)/(cell_counts+1)
            scores=np.repeat((regression[:,:,0]-regression[:,:,1])[:,:,None,None],2,axis=2)
            scores=np.repeat(scores,2,axis=3)
            for arm in [0,1]:
                assignment=treatment if arm==0 else 1-treatment
                for response in [0,1]:
                    scores[:,:,arm,response]+=(1 if arm==0 else -1)*(response-regression[:,:,arm])/assignment
            sums.append(np.sum(evaluation*scores,axis=(1,2,3)))
            squares.append(np.sum(evaluation*scores**2,axis=(1,2,3)))
        estimate=(sums[0]+sums[1])/1000
        variance=(squares[0]+squares[1]-1000*estimate**2)/(1000*999)
        covered=abs(estimate-.1)<=norm.isf(.025)*np.sqrt(variance)
        rows.append(dict(effect_heterogeneity=heterogeneity,repeats=repeats,coverage=float(covered.mean()),
            mcse=float(np.sqrt(covered.mean()*(1-covered.mean())/repeats)),bias=float(np.mean(estimate-.1)),
            rmse=float(np.sqrt(np.mean((estimate-.1)**2))),mean_se=float(np.mean(np.sqrt(variance)))))
    return rows


if __name__=="__main__":
    path=Path(sys.argv[1])
    if path.exists() or not str(path).startswith('/scratch.global/zhan9381/FACE-HD/'):
        raise ValueError("Use a new FACE-HD scratch file")
    result=audit()
    path.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result))
