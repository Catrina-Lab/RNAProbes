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