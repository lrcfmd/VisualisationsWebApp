import matplotlib.pyplot as plt
import numpy as np

class Mcmc:

    def __init__(self,omega,knowns=[],T=1):
        self.omega=omega
        self.knowns=knowns
        self.T=T
    def log_prob(self,positions):
        distances=[]
        for pos in positions:
            for k in self.knowns:
                d=np.linalg.norm(pos-k)
                distances.append(d)
        for n,posa in enumerate(positions[:-1]):
            dist=[]
            for posb in positions[n+1:]:
                d=np.linalg.norm(posa-posb)
                distances.append(d)
        distances=np.array(distances)
        return -self.T*np.sum(1/(distances**2+1e-50))

    def proposal(self,positions,poss_num_swaps,p):
        num_swaps=np.random.choice(poss_num_swaps,p=p)
        positions=positions[np.random.choice(
            len(positions),len(positions)-num_swaps,replace=False)]
        extra=self.omega[
            np.random.choice(len(self.omega),num_swaps,replace=False)]
        return np.concatenate((positions,extra),axis=0)

    def p_acc_MH(self,pos_new,pos_old):
        return min(1,np.exp(self.log_prob(pos_new)-self.log_prob(pos_old)))

    def sample_MH(self,pos_old, poss_num_swaps,p):
        pos_new = self.proposal(pos_old, poss_num_swaps,p)
        # here we determine whether we accept the new state or not:
        # we draw a random number uniformly from [0,1] and compare
        # it with the acceptance probability
        accept = np.random.random() < self.p_acc_MH(pos_new, pos_old)
        if accept:
            return accept, pos_new
        else:
            return accept, pos_old

    def build_MH_chain(self,n_total,n_points,max_swaps=None,plotting=False):
        if max_swaps is None:
            max_swaps=n_points
        init=self.omega[np.random.choice(
            len(self.omega),n_points,replace=False)]
        poss_num_swaps=np.arange(1,max_swaps+1)
        p=poss_num_swaps/sum(poss_num_swaps)
        p=p[::-1]
        n_accepted = 0
        chain = [init]
        score=[self.log_prob(init)]

        for _ in range(n_total):
            accept, state = self.sample_MH(chain[-1],poss_num_swaps,p)
            chain.append(state)
            score.append(self.log_prob(state))
            n_accepted += accept

        acceptance_rate = n_accepted / float(n_total)
        if plotting:
            plt.plot(range(len(score)),score)
            plt.title("Highest score = "+str(round(max(score),2)))
            plt.xlabel("Iteration")
            plt.ylabel("Score")
            plt.show()

        return chain, acceptance_rate, score

'''
model=Model()
model.setup_uncharged(['Li','Al','O'],cube_size=100)
test=Mcmc(model.omega)
score=res[2]
acceptance_rate=res[1]
'''




