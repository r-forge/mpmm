\name{reddeer}
\alias{reddeer}
\docType{data}
\title{
Red Deer dataset.
}
\description{
Red deer engage in aggressive acts to enforce a social hierarchy. The
data consists of T = 1200 observed aggressive actions (e.g. glaring,
kicking, and mounting). Each event had a clear 'winner'
and 'loser', which we code as the sender and recipient, respectively.
}
\usage{data(reddeer)}
\format{
  A data frame with 1200 observations on the following 3 variables.
  \describe{
    \item{\code{winner}}{a factor with levels for each deer that won}
    \item{\code{loser}}{a factor with levels for each deer that lost}
    \item{\code{threat1}}{a factor with levels \code{CT} \code{G1} \code{G2} \code{MA} \code{MG} \code{MO} \code{MR} \code{SA} \code{SW} \code{XX}}
  }
}
\details{
  This is the code sheet as prepared by Michael C. Appleby 
to explain the data in the file REDDEER.txt.  We include this discussion for
reference, and we plan to include the original dataset and our
processing script in later versions of the package.

-----------------------------------------------------

These are all the interactions  recorded from Jan to
March 1978 among stags of known identity in hard horn, i.e.
with 2 normal antlers. Interactions with females, unknown
males and yearlings (whose rudimentary antlers are still
covered with velvet) are excluded. So are those with stags
which had cast one or both antlers: this happened in older
stags first, so older animals appear less and less as March
progresses.

They are not a random sample of interactions from
this bachelor group of stags. They were obtained during a
series of watches on focal individuals, 2 of each age from 2
to 11 years old, as explained in, for example, my 1983 paper
(AB 31, 913). They are thus mostly those involving the focal
individual or members of his party. That same paper discussed
association at some length - which affects availability to be
involved in interactions.

The vast majority of these interactions were
uncontested, so in a sense the fighting terminology of
escalation and winning is inappropriate: hence in the
explanations below I put 'winner' and 'loser' in quotes.
However, we can still ask about who initiates and how
far the 'winner' escalates:

Initiation: there are two ways of getting at this (as
covered in the 1983 paper). First, if one opponent was moving
before the interaction and the other was not, the former must
have initiated. Second, for those interactions which involved
a focal individual I recorded whether an approach to within 5
metres occurred; for the majority which did, the opponent
which approached may be assumed to have initiated the
interaction.
	Escalation: in the 1983 paper I classified interactions into 3 types:

1. Involving only displacement by walking towards the opponent, or very mild threats;
2. With threats such as those listed below as Group 2;
3. With threats or attacks, listed below as Group 3.



% Data columns:
% 1-2     Day
% 4-5     Month
% 7-8     Year
% 10-13   Time
% 15-18   'Winner'
% 20      Antlers of 'winner': all entries in this file should
% 	be H
% 22-25   'Loser' 27 Antlers of 'loser': all entries in this
% 	file should be H
% 29-30   Age of 'winner'
% 32-33   Age of 'loser'
% 35      Prior activity of 'winner':
% 				Moving
% 			Feeding                 
% 			Standing                
% 			Lying           
% 			sParring                
% 			Rubbing                 
% 			Wallowing               
% 			Thrashing               
% 			sNiffing                
% 			X = not recorded
% 37       Prior activity of 'loser'
% 39-40    Vegetation on which 'winner' was standing
% 42-43    Vegetation on which 'loser' was standing
% 45       D = the 'winner' was the focal individual;
% 	S = the 'loser' was the focal individual
% 48-50    Party size (99 = missing value)
% 52-55   Sequence of threats etc.:               
% 			X = none
% 		Group 1:
% 			Ears back or flicked
% 		Group 2:
% 			Glare           
% 		     Antler threats:
% 			Nod             
% 			Wave of antlers                 
% 		     Chin threats:
% 			U = chin up             
% 			High chin threat                
% 			Long reach threat               
% 		  Foot threats:
% 			Foreleg kick (no contact)       
% 		Group 3:
% 			Rush            
% 			Chase           
% 			Prod/push               
% 			K = scissor kick                
% 			Box     
% 		Other:
% 			S = spar                
% 			T = mount               
% 			M = fight
% 57      D = displacement i.e. 'winner' physically displaced
% 		'loser'
% 59      Follow-up:
% 			Feed in place           
% 			C = follow              
% 			Lie in place
% 	Later columns only present for focal individuals
% 63-64   Unimportant
% 66      Approach to within 5 metres:
% 			1 = 'winner'            
% 			2 = 'loser'             
% 			9 = neither
% 67      Subsequent increase of separation to 5 metres or more
% 69-77   Unimportant
}
\references{
M. C. Appleby. Competition in a Red Deer Stag Social Group: Rank, Age and Relatedness of Opponents. Animal Behavior, 31:913-918, 1983.
}
\examples{
data(reddeer)
}
\keyword{datasets}
