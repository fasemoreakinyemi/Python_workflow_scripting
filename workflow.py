import sys
function_names = sys.argv



def hello_world():
	print("hello world")

def good_night():
	print("good night world")

def good_morning():
	print("good morning world")

def main():
	functions = {'hello_world': hello_world,
				'good_night': good_night,
				'good_morning': good_morning}
	if "all" in function_names:
		print("Executing entire workflow")
		for function_name in functions:
			functions[function_name]()
	else:
		print(f"Executing functions {', '.join(sys.argv[1:])}")
		for function_name in sys.argv[1:]:
			functions[function_name]()

main()
