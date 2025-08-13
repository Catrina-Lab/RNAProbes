function createBasicEnum(type, ...Ids){
    let obj = {};
    for(let i = 0; i < Ids.length; i++){
        if(obj[Ids[i]]){ throw new Error('Duplicate IDs are not allowed, including IDs that coincide with default functions')};
        obj[Ids[i]] = Object.freeze({
            toString: () => type + '.' + Ids[i],
            getType: () => type,
            name: Ids[i],
			ordinal: i
        });
    }
	let values = Object.values(obj).sort((a,b)=>a.ordinal - b.ordinal);
    //using defineProperty to prevent counting as enum value
	Object.defineProperty(obj, 'isEnumChild', {
        value: e=>e !== undefined && e === obj[e.name]
    });
	Object.defineProperty(obj, 'getFromOrdinal', {
        value: i=>obj.values[i]
    });
	Object.defineProperty(obj, 'values', {
        value: values
    });
    return Object.freeze(obj);
}

function createElem(type, theClass, parent, extraFuncs = e=> e){
	let elem = document.createElement(type);
	elem.className = theClass;
	extraFuncs(elem);
	parent.appendChild(elem);
	return elem;
}

function wrapIndex(arr, i){
    return arr[((i % arr.length) + arr.length) % arr.length];
}
function getWrapIndex(arr, i){
    return ((i % arr.length) + arr.length) % arr.length;
}

//#region //Validation (assertion) methods
class ValidationError extends Error {
        constructor(message) {
            super(message);
            this.name = 'ValidationError'; // Set a custom name for the error
        }
    }
function validate(boolean, msg){
	if(!boolean){
		throw new ValidationError(msg);
	}
}
//#endregion

//#region //Range and DiscontinousRange
class DiscontinuousRange{
	static isInvalidString(string, minValue, maxValue, forceIncreasing){ //negative so we can return a truthy message
		try{
			DiscontinuousRange.parse(string, minValue, maxValue, forceIncreasing);
			return false;
		} catch (e){
			return e.message || e;
		}
	}
	static parse(rangeStr, minValue, maxValue, forceIncreasing = false){ //forceIncreasing both forces start < stop and n.stop < (b+1).start
		const ranges = rangeStr.replace(/\s+/g, "").split(',');
		return new DiscontinuousRange(minValue, maxValue, forceIncreasing, ...ranges.map(e=>Range.fromString(e)))
	}
    constructor(minValue, maxValue, forceIncreasing, ...items){
        if(items.length == 0){
			throw new Error("You must include a list of ranges to use in DiscontinousRange. Use DiscontinousRange.parse(str) to get one from a string");
		}
		const bound = new Bound(minValue, maxValue);
		validate(items.every(range=>bound.containsRange(range)), `One of the ranges is out of bounds. All must be between ${bound.start} inclusive and ${bound.stop} exclusive`);
        if(forceIncreasing) validate(items.every((range, i)=>(range.start < range.stop) && (i == items.length - 1 || range.stop <= items[i+1].start)), `The given ranges must be increasing (both themselves and relative to each other)`);
		this.ranges = items;
	}
    
    static template(inputString = true, minValue, maxValue){
        if(inputString) return x=> new DiscontinuousRange.parse(x, minValue, maxValue)
        else return (...x) => new DiscontinuousRange(minValue, maxValue, ...x)
	}
	*generator(){
		for(let range of this.ranges){
			for(let i of range){
				yield i;
			}
		}
	}
	[Symbol.iterator](){
		return this.generator();
	}
}
class Range{
	static fromString(string, delim = ":"){
		const portions = string.split(delim);
		const endI = portions.length - 1;
		if((portions[0] !== "" && isNaN(portions[0])) || (portions[endI] !== "" && isNaN(portions[endI]))){
			throw new Error(`Invalid Range String. It should use a "${delim}" as a separator and consist of only numbers (or nothing)`);
		}
		const start = portions[0] === "" ? -Infinity : parseInt(portions[0]);
		const stop = portions[endI] === "" ? Infinity : parseInt(portions[endI]) + 1;
		return new Range(start, stop);
	}
	constructor(start, stop, step = 1){
		this.start = start ?? -Infinity;
		this.stop = stop ?? Infinity;
		this.step = step;
		if(Math.abs(this.start) === Infinity && Math.abs(this.step) !== 1){
			throw new TypeError("Can't have an infinite start without a step magnitude of 1");
		}
		if(this.step == 0){
			throw new TypeError("Can't have an Range step of 0");
		}
	}
	trueEnd(){
		if(this.stop == Infinity){
			return Infinity;
		}
		if(this.step > 0 ? this.start >= this.stop : this.stop >= this.start){
			return NaN;
		}
		const inclusive = this.stop - Math.sign(this.step);
		const offset = (this.start === Infinity ? 0 : ((inclusive - this.start) % this.step));
		return inclusive - offset
	}
	containsRange(range){
		throw new Error("Not implemented yet");
		thisTrueEnd = this.trueEnd();
		testTrueEnd = range.trueEnd();

		if(testTrueEnd === NaN) return true;
		if(thisTrueEnd === NaN) return false;

		
	}
	*generator(){
		if(Math.abs(this.start) === Infinity || Math.abs(this.stop) === Infinity){
			throw new Error("Can't get a generator for an infinite range");
		}
		for(let i = this.start; this.step > 0 ? i < this.stop : i > this.stop; i+=this.step){
			yield i;
		}
	}
	[Symbol.iterator](){
		return this.generator();
	}
}
class Bound extends Range{
	constructor(start, stop){
		super(start, stop, 1);
		if(start > stop){
			throw new TypeError("Bound can't have a start > stop");
		}
	}
	containsValue(val){
		return isInBound(this.start, this.stop, val);
	}
	containsRange(rangeObject){
		validate(rangeObject instanceof Range, "Input must be of type Range")
		let rStart = rangeObject.start, rStop = rangeObject.stop - 1 //actual range values
		if(rangeObject.step < 1){
			 rStart = rangeObject.stop+1, rStop = rangeObject.start
		}
		return isNotBelowMin(this.start, rStart) && isBelowMax(this.stop, rStop)
	}
}
function isNotBelowMin(min, ...values){
	return min === undefined || values.every(value => value != undefined && value >= min);
}
function isBelowMax(max, ...values){
	return max === undefined || values.every(value => value != undefined && value < max);
}

function isInBound(min, max, ...values){
    return isNotBelowMin(min, ...values) && isBelowMax(max, ...values);
}

//#endregion