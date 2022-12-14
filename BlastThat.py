import random
import numpy as np
import matplotlib.pyplot as plt


class Blaster(object):
    def __init__(self, T):
        self.T_data = T
        self.S_data = ""
        self.scores = None
        self.score_matrix = None
        self.s_solution = ""
        self.t_solution = ""
        self.alignementcoordo = [0,0]
        
    def print2(self):
        print(self.T_data)
        
    def print3(self, P):
        print(self.T_data)
        print(P)
    
    def blast_first_pass(self, S_sq):
        self.scores = np.zeros(len(self.T_data))
        indexes = np.array(range(len(self.T_data)))
        
        S_mid = int(np.floor(len(S_sq) / 2))
        small_S_sq = S_sq[S_mid-1: S_mid+2]
        print(small_S_sq)
        self.blast_this(small_S_sq,indexes, 1)
            #for j in range(-S_mid, S_mid):
            
    def blast_second_pass(self, S_sq):
        maxscore = np.amax(self.scores)
        indexes = np.argwhere(self.scores == maxscore)
        indexes = np.squeeze(indexes)
        self.scores = np.zeros(len(self.T_data))
        self.scores = self.scores - len(S_sq)
        self.blast_this(S_sq, indexes, len(S_sq)/2)
        
    
    def blast_this(self,S_sq, index_to_look_T, each_side_step, should_penalise_array = True):
        if (self.scores is None):
            self.scores = np.zeros(len(self.T_data))
            
        for i in index_to_look_T:
            S_mid = int(np.floor(len(S_sq) / 2))
            hamming_dist = int(self.T_data[i] != S_sq[S_mid])
            dstep = 1
            while hamming_dist < 10:
                di = i - dstep
                dj = S_mid - dstep
                if (di < 0 or dj < 0):
                    hamming_dist += 1
                else:
                    hamming_dist += int(self.T_data[di] != S_sq[dj])
                    
                di = i + dstep
                dj = S_mid + dstep
                if (di >= len(self.T_data) or dj >= len(S_sq)):
                    if should_penalise_array:
                        hamming_dist += 1
                else:
                    hamming_dist += int(self.T_data[di] != S_sq[dj])
                
                dstep +=1
                if (dstep > each_side_step):
                    break
            
            self.scores[i] = -hamming_dist
        return self.scores
            
    def full_blast(self, S_sq):
        blast_scores = np.zeros([len(S_sq), len(self.T_data)])
        self.S_data = S_sq
        s_len = len(S_sq)
        
        for j in range(s_len):
            for i in range(len(self.T_data)):
                score = 0
                if self.T_data[i] == S_sq[j]:
                    score = 1
                
                min_i = i
                min_j = i
                max_i = j
                max_j = j
                dstep = 1
                while score > 0 and dstep < 3:
                    di = i - dstep
                    dj = j - dstep
                    if (di < 0 or dj < 0 or self.T_data[di] != S_sq[dj]):
                        score -= 1
                    else:
                        score +=1
                        min_i = min(di, min_i)
                        min_j = min(dj, min_j)
                    
                    di = i + dstep
                    dj = j + dstep
                    if (di >= len(self.T_data) or dj >= len(S_sq) or self.T_data[di] != S_sq[dj]):
                        score -= 1
                    else:
                        score +=1
                        max_i = max(di, max_i)
                        max_j = max(dj, max_j)
                        
                    dstep +=1
                
                di = min_i
                dj = min_j
                while di <= max_i and dj <= max_j:
                    if di < 0 or dj < 0:
                        pass
                    if di >= len(self.T_data) or dj >= len(S_sq):
                        pass
                    blast_scores[dj, di] = max(score, blast_scores[dj, di])
                    di+=1
                    dj+=1
                    
        self.score_matrix = blast_scores
        
    def matrixscore(self, i, j) -> float:
        if i < 0 or j < 0:
            return 0.0
        if i >= len(self.score_matrix[0]) or j >= len(self.score_matrix):
            return 0.0
        
        return self.score_matrix[j,i]
    
    def searchTopLeft(self, i, j, reach):
        cur = self.matrixscore(i,j)
        
        if reach <= 0:
            return cur, i ,j
        
        top, ti, tj = self.searchTopLeft(i, j-1, reach-1)
        left, li, lj = self.searchTopLeft(i-1, j, reach-1)
        if (top > cur and top > left):
            self.score_matrix[j, i] = top
            return top, ti, tj
        elif(left > cur and left > top):
            self.score_matrix[j, i] = left
            return left, li, lj
        else:
            return cur, i, j
        
    def searchBotRight(self, i, j, reach):
        cur = self.matrixscore(i,j)
        
        if reach <= 0:
            return cur, i ,j
        
        bot, ti, tj = self.searchBotRight(i, j+1, reach-1)
        right, li, lj = self.searchBotRight(i+1, j, reach-1)
        if (bot > cur and bot > right):
            self.score_matrix[j, i] = bot
            return bot, ti, tj
        elif(right > cur and right > bot):
            self.score_matrix[j, i] = right
            return right, li, lj
        else:
            return cur, i, j

    def extend_lines(self, treshold, search_reach):
        mediane_line = int(len(self.S_data)/2)
        longest = 0
        endline = [0, 0]
        for i in range(len(self.T_data)):
            cur = self.matrixscore(i, mediane_line)
            if cur > treshold:
                line_len = 0
                # vers le haut
                di = i
                dj = mediane_line
                last = [dj, di]
                while cur > treshold:

                    line_len += 1
                    di -= 1
                    dj -= 1
                    cur = self.matrixscore(di, dj)
                    if cur <= treshold:
                        cur, di, dj = self.searchTopLeft(di, dj, search_reach)

                
                # vers le bas
                di = i+1
                dj = mediane_line+1
                cur = self.matrixscore(di, dj)
                while cur > treshold:
                    line_len += 1
                    last = [dj, di]
                    di += 1
                    dj += 1
                    cur = self.matrixscore(di, dj)
                    if cur <= treshold:
                        cur, di, dj = self.searchBotRight(di, dj, search_reach)

                
                if line_len > longest:
                    longest = line_len
                    endline = last

        self.alignementcoordo = endline
        
    def trouver_sequence_allignee(self, treshold, search_reach):
        S = "   "
        T = ""

        line_len = 0
        j = self.alignementcoordo[0]
        i = self.alignementcoordo[1]
        cur = self.matrixscore(i,j)
        if cur <= treshold:
            print(cur)
            print(self.alignementcoordo)
        
        while cur > treshold:

            T += self.T_data[i]
            S += self.S_data[j]
            line_len += 1
            i -= 1
            j -= 1
            cur = self.matrixscore(i, j)
            if cur <= treshold:
                old_i = i
                old_j = j
                cur, i, j = self.searchTopLeft(i, j, search_reach)
                if (old_i != i and old_j == j):
                    T += "-"
                elif (old_i == i and old_j == j):
                    S += "-"
                elif (old_i != i and old_j != j):
                    T += self.T_data[i]
                    S += self.S_data[j]
        
        sol_T = self.T_data[self.alignementcoordo[1]-4:self.alignementcoordo[1]-1]
        sol_T += T
        sol_T += self.T_data[i+1:i+4]
        print('T', sol_T)
        print('S', S)
    
            
                
    def print_blast_matrix(self):
        plt.figure()
        plt.subplots(figsize=(40,40))
        plt.imshow(self.score_matrix, cmap='autumn',interpolation='nearest')

    def print_scores(self):
        xdata = np.arange(1, len(self.scores) + 1)
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.ylabel('score')
        plt.plot(xdata, self.scores, label='Sequence Score')
        plt.legend()
