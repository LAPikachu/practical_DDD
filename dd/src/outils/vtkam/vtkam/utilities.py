def get_string(prompt):
    while True:
        try:
            valtmp = raw_input(prompt)
            if (valtmp =='r') :
              value=-1
              break
            value =  map(str,valtmp.split(' '))
            if len(value) != 1:
              print("Sorry, 1 str value is expected.")
              continue
            else:
              break
        except ValueError:
            print("Sorry, I didn't understand that.")
            continue
        except KeyboardInterrupt: quit()

    return value

def get_int(prompt):
    while True:
        try:
            valtmp = raw_input(prompt)
            if (valtmp =='r') :
              value=-1
              break
            value = int(valtmp)
        except ValueError:
            print("Sorry, I didn't understand that.")
            continue
        except KeyboardInterrupt: quit()
        except : break

    return value

def get_non_negative_int(prompt):
    while True:
        try:
            valtmp = raw_input(prompt)
            if (valtmp =='r') :
              value=-1
              break
            value = int(valtmp)
            if value < 0:
		  				print("Sorry, your response must not be negative.")
			  			continue
            else:
              break
        except ValueError:
            print("Sorry, I didn't understand that.")
            continue
        except KeyboardInterrupt: quit()
        except : break

    return value

def get_int3(prompt):
    while True:
        try:
            valtmp = raw_input(prompt)
            if (valtmp =='r') :
              value=-1
              break
            values = map(int,valtmp.split(' '))
            if len(values) != 3:
              print("Sorry, 3 integer values are expected.")
              continue
            else:
              break
        except ValueError:
            print("Sorry, I didn't understand that.")
            continue
        except KeyboardInterrupt: quit()
        except : break

    return values

def get_non_negative_float(prompt):
    while True:
        try:
            valtmp = raw_input(prompt)
            if (valtmp =='r') :
              value=-1
              break
            value = float(valtmp)
            if value < 0:
              print("Sorry, your response must not be negative.")
              continue
            else:
              break
        except ValueError:
            print("Sorry, I didn't understand that.")
            continue
        except KeyboardInterrupt: quit()
        except : break

    return value

def get_float(prompt):
    while True:
        try:
            valtmp = raw_input(prompt)
            if (valtmp =='r') :
              value=-1
              break
            value = float(valtmp)
        except ValueError:
            print("Sorry, I didn't understand that.")
            continue
        except KeyboardInterrupt: quit()
        except: break

    return value

def get_float3(prompt):
    while True:
        try:
            valtmp = raw_input(prompt)
            if (valtmp =='r') :
              value=-1
              break
            values = map(float,valtmp.split(' '))
            if len(values) != 3:
              print("Sorry, 3 integer values are expected.")
              continue
            else:
              break
        except ValueError:
            print("Sorry, I didn't understand that.")
            continue
        except KeyboardInterrupt: quit()
        except : break

    return values