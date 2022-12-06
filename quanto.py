from quantum import *

rps = qprogram(2)

rps.addgates(0, [HGATE(), RGATE(pi/3.3)])
rps.addgates(1, [IGATE(), CNOTGATE(), IGATE()])

rps.compile(showcompilationresult = False)

beats = {
    'r' : 's',
    's' : 'p',
    'p' : 'r'
}

emojis = {
    'r' : '💎',
    'p' : '📃',
    's' : '✂'
}

while True:
    userinput = input("Enter your guess\n[(q)uit or (e)xit to stop]\n(R)ock, (P)aper or (S)cissors : ").strip().lower()

    if userinput in ['exit', 'quit', 'q', 'e', ''] : break
    userinput = userinput[0]

    if userinput in ['r', 'p', 's']:
        compchoice = ['r', 'p', 'p', 's'][rps.measure().index(1)]

        print(f"\nYou -> {emojis[userinput]} {emojis[compchoice]} <- Quanto")
        
        if beats[userinput] == compchoice: print("You Win!")
        elif userinput == compchoice: print("It's a draw!")
        else : print("You lose!")

    else : print("Invalid choice!\n")