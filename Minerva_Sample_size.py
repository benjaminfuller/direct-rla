from math import ceil, sqrt
from scipy.stats import norm

@staticmethod
def round_size_approx(margin, alpha, quant):
        """
        Returns approximate round size for small margins
        :param margin: margin of victory (float in [0, 1])
        :param alpha: risk limit
        :param quant: desired probability of stopping in the next round
        :return: the next round size computed under a normal approximation to the binomial
        """
        z_a = norm.isf(quant)
        z_b = norm.isf(alpha * quant)
        p = (1 + margin) / 2
        return ceil(((z_a * sqrt(p * (1 - p)) - .5 * z_b) / (p - .5)) ** 2)



if __name__ == "__main__":
        quant = .9
        alpha = .05

        margins = [0.03, 0.025, .02, .019, .018, .017, .016, .015, .014, .013, .012, .011, .01, .009, .008,.0075, .007, .006, .005, .004, .003, .002]

        
        size_dict = {}
        for margin in margins:
                round_size =round_size_approx(margin, alpha, quant)
                size_dict[margin] = round_size
        
        print(size_dict)
